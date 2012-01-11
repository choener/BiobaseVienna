{-# LANGUAGE DoAndIfThenElse #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE NoMonomorphismRestriction #-}

-- | Import ViennaRNA parameter files (*par files). Imports Turner2004 files.

module Biobase.Vienna.ImportPar where

import Data.Iteratee as I
import Data.Iteratee.Char as I
import Data.Iteratee.IO as I
import Data.ByteString.Char8 as BS
import Prelude as P
import Data.Array.Repa.Index
import Data.Array.Repa.Shape
import qualified Data.Map as M
import Control.Monad.Identity

import Data.PrimitiveArray
import Data.PrimitiveArray.Unboxed
import Biobase.Secondary.Vienna
import Biobase.Primary

import Biobase.Vienna
import Biobase.Vienna.Import



-- | split the input into different blocks. Each block has a name (prefixed #
-- in the file) and some data.


data BL
  = Block {fromBlock :: [Int]}
  | Lookup {fromLookup :: [(ByteString,Int,Int)]}
  deriving (Show)

lookupStructure =
  [ "# Triloops"
  , "# Tetraloops"
  , "# Hexaloops"
  ]

iterBlocks = enumLinesBS ><> I.filter (not . BS.null) ><> convStream f where
  f = do
    h <- I.head -- the # prefix
    if h `P.elem` lookupStructure
    then do
      xs' <- I.takeWhile ((/='#') . BS.head)
      let xs = P.map (\[w,x,y] -> (w,getInt x, getInt y)) . P.map BS.words $ xs'
      return [(h, Lookup xs)]
    else do
      xs' <- I.takeWhile ((/='#') . BS.head)
      let xs = P.map getInt . P.concat . P.map BS.words . P.map (BS.takeWhile (/='/')) $ xs'
      return [(h, Block xs)]

getInt :: ByteString -> Int
getInt s
  | s == "INF" = 999999
  | otherwise  = read . BS.unpack $ s

fromFile :: FilePath -> IO (Vienna2004,Vienna2004)
fromFile fp = do
  i <- enumFile 8192 fp (joinI $ iterBlocks stream2list)
  bs <- run i
  return $ makeStructures bs

fromByteString :: ByteString -> (Vienna2004,Vienna2004)
fromByteString s = runIdentity $ do
  i <- enumPure1Chunk s (joinI $ iterBlocks stream2list)
  bs <- run i
  return $ makeStructures bs

makeStructures bs =
  let vEner = Vienna2004
        { stack = blockAssocs minPP maxPP ppKeys   $ lookup "# stack"   bs
        , dangle3 = blockAssocs minPB maxPB pbKeys $ lookup "# dangle3" bs
        , dangle5 = blockAssocs minPB maxPB pbKeys $ lookup "# dangle5" bs
        , hairpinL = blockAssocs (Z:.0) (Z:.30) thirty $ lookup "# hairpin" bs
        , hairpinMM = blockAssocs minPBB maxPBB pbbKeys $ lookup "# mismatch_hairpin" bs
        , hairpinLookup = fst $ allLookups bs
        , hairpinGGG = 999999
        , hairpinCslope = 999999
        , hairpinCintercept = 999999
        , hairpinC3 = 999999
        , bulgeL = blockAssocs (Z:.0) (Z:.30) thirty $ lookup "# bulge" bs
        , bulgeSingleC = 999999
        , iloop1x1 = blockAssocs minPPBB maxPPBB ppbbKeys $ lookup "# int11" bs
        , iloop2x1 = blockAssocs minPPBBB maxPPBBB ppbbbKeys $ lookup "# int21" bs
        , iloop2x2 = blockAssocs minPPBBBB maxPPBBBB ppbbbbKeys $ lookup "# int22" bs
        , iloopMM = blockAssocs minPBB maxPBB pbbKeys $ lookup "# mismatch_interior" bs
        , iloop2x3MM = blockAssocs minPBB maxPBB pbbKeys $ lookup "# mismatch_interior_23" bs
        , iloop1xnMM = blockAssocs minPBB maxPBB pbbKeys $ lookup "# mismatch_interior_1n" bs
        , iloopL = blockAssocs (Z:.0) (Z:.30) thirty $ lookup "# interior" bs
        , multiMM = blockAssocs minPBB maxPBB pbbKeys $ lookup "# miscmatch_multi" bs
        , ninio = 999999
        , maxNinio = 999999
        , multiOffset = 999999
        , multiNuc = 999999
        , multiHelix = 999999
        , multiAsym = 999999
        , multiStrain = 999999
        , extMM = blockAssocs minPBB maxPBB pbbKeys $ lookup "# mismatch_exterior" bs
        , coaxial = fromAssocs minPP maxPP 999999 []
        , coaxStack = fromAssocs minPBB maxPBB 999999 []
        , tStackCoax = fromAssocs minPBB maxPBB 999999 []
        , largeLoop = 999999
        , termAU = 999999
        , intermolecularInit = 999999
        }
      vEnth = Vienna2004
        {}
  in (vEner,vEnth)

thirty = P.map (Z:.) [0..30]
pbKeys = [ (Z:.x:.y) | x<-cgnsP, y<-nacgu ]
ppKeys = [ (Z:.x:.y) | x<-cgnsP, y<-cgnsP ]
pbbKeys = [ (Z:.a:.b:.c) | a<-cgnsP, b<-nacgu, c<-nacgu ]
ppbbKeys = [ pp:.a:.b | pp<-ppKeys, a<-nacgu, b<-nacgu ]
ppbbbKeys = [ ppbb:.a | ppbb<-ppbbKeys, a<-nacgu ]
ppbbbbKeys = [ Z:.a:.b:.c:.d:.e:.f | a<-cguaP, b<-cguaP, c<-acgu, d<-acgu, e<-acgu, f<-acgu ]
blockAssocs minKey maxKey keys (Just (Block xs)) = fromAssocs minKey maxKey 999999 $ P.zip keys xs
allLookups bs = (M.fromList $ P.map (\(a,b,c) -> (a,b)) ls, M.fromList $ P.map (\(a,b,c) -> (a,c)) ls) where
  ls = P.map (\(a,b,c) -> (mkPrimary a,b,c)) $ P.concatMap (fromLookup . snd) $ P.filter (isL . snd) bs
  isL (Lookup _) = True
  isL _          = False
