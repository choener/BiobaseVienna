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

import Data.PrimitiveArray
import Data.PrimitiveArray.Unboxed
import Biobase.Secondary.Vienna
import Biobase.Primary

import Biobase.Vienna
import Biobase.Vienna.Import



fromFile :: FilePath -> IO (Vienna2004,Vienna2004)
fromFile fp = do
  return undefined

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

--iterBlocks :: Iteratee ByteString m [(ByteString,Either Block Lookup)]
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

test = do
  i <- enumFile 8192 "parfiles/rna_turner2004.par" (joinI $ iterBlocks stream2list)
  bs <- run i
  let vEner = Vienna2004
        { stack = blockAssocs minPP maxPP ppKeys   $ lookup "# stack"   bs
        , dangle3 = blockAssocs minPB maxPB pbKeys $ lookup "# dangle3" bs
        , dangle5 = blockAssocs minPB maxPB pbKeys $ lookup "# dangle5" bs
        , hairpinL = blockAssocs (Z:.0) (Z:.30) thirty $ lookup "# hairpin" bs
        , hairpinMM = blockAssocs minPBB maxPBB pbbKeys $ lookup "# mismatch_hairpin" bs
        , hairpinLookup :: Map [Nuc] Int
        , hairpinGGG :: Int
        , hairpinCslope :: Int
        , hairpinCintercept :: Int
        , hairpinC3 :: Int
        }
  print . assocs $ dangle5 vEner
  return vEner

thirty = P.map (Z:.) [0..30]
pbKeys = [ (Z:.x:.y) | x<-cgnsP, y<-nacgu ]
ppKeys = [ (Z:.x:.y) | x<-cgnsP, y<-cgnsP ]
pbbKeys = [ (Z:.a:.b:.c) | a<-cgnsP, b<-nacgu, c<-nacgu ]
blockAssocs minKey maxKey keys (Just (Block xs)) = fromAssocs minKey maxKey 999999 $ P.zip keys xs
