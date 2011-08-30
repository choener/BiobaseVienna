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
import qualified Data.ByteString.Char8 as BS

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
  , hairpinLookup = M.mapKeys (read . BS.unpack) . M.map deka $ T.hairpinLookup turner
  , hairpinGGG = deka $ T.hairpinGGG turner
  , hairpinCslope = deka $ T.hairpinCslope turner
  , hairpinCintercept = deka $ T.hairpinCintercept turner
  , hairpinC3 = deka $ T.hairpinC3 turner
  , bulgeL = convert 0 30 $ T.bulgeL turner
  , bulgeSingleC = deka $ T.bulgeSingleC turner
  , iloop1x1 = convert minPPBB maxPPBB $ T.iloop1x1 turner
  , iloop2x1 = convert minPPBBB maxPPBBB $ T.iloop2x1 turner
  , iloop2x2 = convert minPPBBBB maxPPBBBB $ T.iloop2x2 turner
  , iloopMM = convert minPBB maxPBB $ T.iloopMM turner
  , iloop2x3MM = convert minPBB maxPBB $ T.iloop2x3MM turner
  , iloop1xnMM = convert minPBB maxPBB $ T.iloop1xnMM turner
  , iloopL = convert 0 30 $ T.iloopL turner
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

instance (IdxConvert al ar, IdxConvert bl br, IdxConvert cl cr, IdxConvert dl dr) => IdxConvert (al,bl,cl,dl) (ar,br,cr,dr) where
  idxConvert (a,b,c,d) = (idxConvert a, idxConvert b, idxConvert c, idxConvert d)

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
minPPBB = (minP,minP,(nN,nN))
maxPPBB = (maxP,maxP,(nU,nU))
minPPBBB = (minP,minP,(nN,nN,nN))
maxPPBBB = (maxP,maxP,(nU,nU,nU))
minPPBBBB = (minP,minP,(nN,nN,nN,nN))
maxPPBBBB = (maxP,maxP,(nU,nU,nU,nU))

