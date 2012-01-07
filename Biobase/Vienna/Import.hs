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
import Data.Array.Repa.Index
import Data.Array.Repa.Shape
import qualified Data.ByteString.Char8 as BS
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
  , hairpinL = convert (Z:.0) (Z:.30) $ T.hairpinL turner
  , hairpinMM = convert minPBB maxPBB $ T.hairpinMM turner
  , hairpinLookup = M.mapKeys (read . BS.unpack) . M.map deka $ T.hairpinLookup turner
  , hairpinGGG = deka $ T.hairpinGGG turner
  , hairpinCslope = deka $ T.hairpinCslope turner
  , hairpinCintercept = deka $ T.hairpinCintercept turner
  , hairpinC3 = deka $ T.hairpinC3 turner
  , bulgeL = convert (Z:.0) (Z:.30) $ T.bulgeL turner
  , bulgeSingleC = deka $ T.bulgeSingleC turner
  , iloop1x1 = convert minPPBB maxPPBB $ T.iloop1x1 turner
  , iloop2x1 = convert minPPBBB maxPPBBB $ T.iloop2x1 turner
  , iloop2x2 = convert minPPBBBB maxPPBBBB $ T.iloop2x2 turner
  , iloopMM = convert minPBB maxPBB $ T.iloopMM turner
  , iloop2x3MM = convert minPBB maxPBB $ T.iloop2x3MM turner
  , iloop1xnMM = convert minPBB maxPBB $ T.iloop1xnMM turner
  , iloopL = convert (Z:.0) (Z:.30) $ T.iloopL turner
  , multiMM = convert minPBB maxPBB $ T.multiMM turner
  , ninio = deka $ T.ninio turner
  , maxNinio = deka $ T.maxNinio turner
  , multiOffset = deka $ T.multiOffset turner
  , multiNuc = deka $ T.multiNuc turner
  , multiHelix = deka $ T.multiHelix turner
  , multiAsym = deka $ T.multiAsym turner
  , multiStrain = deka $ T.multiStrain turner
  , extMM = convert minPBB maxPBB $ T.extMM turner
  , coaxial = convert minPP maxPP $ T.coaxial turner
  , coaxStack = convert minPBB maxPBB $ T.coaxStack turner
  , tStackCoax = convert minPBB maxPBB $ T.tStackCoax turner
  , largeLoop = deka $ T.largeLoop turner
  , termAU = deka $ T.termAU turner
  , intermolecularInit = deka $ T.intermolecularInit turner
  } where convert mn mx = fromAssocs mn mx 999999 . map (idxConvert *** deka) . assocs

-- conversion of indices
--
-- (roll your eyes)

class IdxConvert a b where
  idxConvert :: a -> b

instance IdxConvert DIM1 DIM1 where
  idxConvert = id

instance IdxConvert T.PNN PNN where
  idxConvert (Z:.p1:.p2:.n1:.n2) = Z:. mkViennaPair (p1,p2) :.n1:.n2

instance IdxConvert T.PN PN where
  idxConvert (Z:.p1:.p2:.n1) = Z:. mkViennaPair (p1,p2) :.n1

instance IdxConvert T.PPNNNN PPNNNN where
  idxConvert (Z:.p11:.p12:.p21:.p22:.n1:.n2:.n3:.n4) = Z:. mkViennaPair (p11,p12) :. mkViennaPair (p21,p22) :.n1:.n2:.n3:.n4

instance IdxConvert T.PPNNN PPNNN where
  idxConvert (Z:.p11:.p12:.p21:.p22:.n1:.n2:.n3) = Z:. mkViennaPair (p11,p12) :. mkViennaPair (p21,p22) :.n1:.n2:.n3

instance IdxConvert T.PPNN PPNN where
  idxConvert (Z:.p11:.p12:.p21:.p22:.n1:.n2) = Z:. mkViennaPair (p11,p12) :. mkViennaPair (p21,p22) :.n1:.n2

instance IdxConvert T.PP PP where
  idxConvert (Z:.p11:.p12:.p21:.p22) = Z:. mkViennaPair (p11,p12) :. mkViennaPair (p21,p22)



-- | Transform energies to the vienna Int-based variant
--
-- (which is round (e*100)).

deka = round . (*100)

minP = Z:.vpNP -- minBound
maxP = Z:.vpNS -- maxBound
minPB = minP:.nN -- (minP,nN)
maxPB = maxP:.nU -- (maxP,nU)
minPP = minP:.vpNP -- (minP,minP)
maxPP = maxP:.vpNS -- (maxP,maxP)
minPBB = minP:.nN:.nN -- (minP,nN,nN)
maxPBB = maxP:.nU:.nU -- (maxP,nU,nU)
minPPBB = minPP:.nN:.nN -- (minP,minP,(nN,nN))
maxPPBB = maxPP:.nU:.nU -- (maxP,maxP,(nU,nU))
minPPBBB = minPPBB:.nN -- (minP,minP,(nN,nN,nN))
maxPPBBB = maxPPBB:.nU -- (maxP,maxP,(nU,nU,nU))
minPPBBBB = minPPBBB:.nN -- (minP,minP,(nN,nN,nN,nN))
maxPPBBBB = maxPPBBB:.nU -- (maxP,maxP,(nU,nU,nU,nU))

