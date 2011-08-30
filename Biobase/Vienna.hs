
-- | The Vienna 2004 model is a mirror of the Turner data structure. The
-- difference is that The ViennaRNA package uses an explicit annotation for RNA
-- basepairs.
--
-- Furthermore, all energies are stored in parts of 1/100, using Ints instead
-- of Doubles.

module Biobase.Vienna where

import Data.ByteString (ByteString)
import Data.Map (Map)

import Data.PrimitiveArray
import Biobase.Primary
import Biobase.Secondary.Vienna

data Vienna2004 = Vienna2004
  { stack :: PrimArray (ViennaPair,ViennaPair) Int
  , dangle3 :: PrimArray PN Int
  , dangle5 :: PrimArray PN Int
  , hairpinL :: PrimArray Int Int
  , hairpinMM :: PrimArray PNN Int
  , hairpinLookup :: Map [Nuc] Int
  , hairpinGGG :: Int
  , hairpinCslope :: Int
  , hairpinCintercept :: Int
  , hairpinC3 :: Int
  , bulgeL :: PrimArray Int Int
  , bulgeSingleC :: Int
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

type PNN = (ViennaPair,Nuc,Nuc)
type PN  = (ViennaPair,Nuc)
