{-# LANGUAGE TypeSynonymInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FlexibleInstances #-}

-- | There are two import scenarios: (i) importing from a ViennaRNA parameter
-- file (version 2004 only) and (ii) importing from a Turner2004 data
-- structure.
--
-- TODO (i) is missing right now

module Biobase.Vienna.Import where

import Control.Arrow
import qualified Data.Map as M

import Biobase.Primary
import Biobase.Secondary
import Biobase.Secondary.Vienna
import Data.PrimitiveArray
import qualified Biobase.Turner as T

import Biobase.Vienna



-- * Transforming a Turner2004 data structure into a Vienna2004 data structure.

-- | From a 'Turner2004' data structure via lists of key/value pairs.

fromTurner2004 :: T.Turner2004 -> Vienna2004
fromTurner2004 turner = Vienna2004
  { stack = convert minPP maxPP $ T.stack turner
  , dangle3 = convert minPB maxPB $ T.dangle3 turner
  , dangle5 = convert minPB maxPB $ T.dangle5 turner
  , hairpinL = convert 0 30 $ T.hairpinL turner
  , hairpinMM = convert minPBB maxPBB $ T.hairpinMM turner
  , hairpinLookup = M.map deka $ T.hairpinLookup turner
  , hairpinGGG = deka $ T.hairpinGGG turner
  , hairpinCslope = deka $ T.hairpinCslope turner
  , hairpinCintercept = deka $ T.hairpinCintercept turner
  , hairpinC3 = deka $ T.hairpinC3 turner
  , bulgeL = convert 0 30 $ T.bulgeL turner
  , bulgeSingleC = deka $ T.bulgeSingleC turner
  } where convert mn mx = fromAssocs mn mx 999999 . map (idxConvert *** deka) . assocs

-- conversion of indices
--
-- (roll your eyes)

class IdxConvert a b where
  idxConvert :: a -> b

instance IdxConvert Pair ViennaPair where
  idxConvert = mkViennaPair

instance IdxConvert Nuc Nuc where
  idxConvert = id

instance IdxConvert Int Int where
  idxConvert = id

instance (IdxConvert a c, IdxConvert b d) => IdxConvert (a,b) (c,d) where
  idxConvert = idxConvert *** idxConvert

instance (IdxConvert a d, IdxConvert b e, IdxConvert c f) => IdxConvert (a,b,c) (d,e,f) where
  idxConvert (a,b,c) = (idxConvert a, idxConvert b, idxConvert c)

-- | Transform energies to the vienna Int-based variant
--
-- (which is round (e*100)).

deka = round . (*100)

minP = minBound
maxP = maxBound
minPB = (minP,nN)
maxPB = (maxP,nU)
minPP = (minP,minP)
maxPP = (maxP,maxP)
minPBB = (minP,nN,nN)
maxPBB = (maxP,nU,nU)

{-
  , iloop1x1 :: PrimArray (ViennaPair,ViennaPair,(Nuc,Nuc)) Int
  , iloop2x1 :: PrimArray (ViennaPair,ViennaPair,(Nuc,Nuc,Nuc)) Int
  , iloop2x2 :: PrimArray (ViennaPair,ViennaPair,(Nuc,Nuc,Nuc,Nuc)) Int
  , iloopMM :: PrimArray PNN Int
  , iloop2x3MM :: PrimArray PNN Int
  , iloop1xnMM :: PrimArray PNN Int
  , iloopL :: PrimArray Int Int
  , multiMM :: PrimArray PNN Int
  , ninio :: Int
  , maxNinio :: Int
  , multiOffset :: Int
  , multiNuc :: Int
  , multiHelix :: Int
  , multiAsym :: Int
  , multiStrain :: Int
  , extMM :: PrimArray PNN Int
  , coaxial :: PrimArray (ViennaPair,ViennaPair) Int -- no intervening unpaired nucleotides
  , coaxStack :: PrimArray PNN Int
  , tStackCoax :: PrimArray PNN Int
  , largeLoop :: Int
  , termAU :: Int
  , intermolecularInit :: Int
  } -- deriving (Read,Show)
-}
