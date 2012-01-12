{-# LANGUAGE TypeOperators #-}

-- | The Vienna 2004 model is a mirror of the Turner data structure. The
-- difference is that The ViennaRNA package uses an explicit annotation for RNA
-- basepairs.
--
-- Furthermore, all energies are stored in parts of 1/100, using Ints instead
-- of Doubles.

module Biobase.Vienna where

import Data.Array.Repa.Index
import Data.Array.Repa.Shape
import Data.ByteString (ByteString)
import Data.Map (Map)

import Data.PrimitiveArray
import Biobase.Primary
import Biobase.Secondary.Vienna

type P = Z:.ViennaPair
type PN = P:.Nuc
type PNN = PN:.Nuc
type PP = Z:.ViennaPair:.ViennaPair
type PPNN = PP:.Nuc:.Nuc
type PPNNN = PPNN:.Nuc
type PPNNNN = PPNNN:.Nuc

data Vienna2004 = Vienna2004
  { stack :: PrimArray PP Int
  , dangle3 :: PrimArray PN Int
  , dangle5 :: PrimArray PN Int
  , hairpinL :: PrimArray DIM1 Int
  , hairpinMM :: PrimArray PNN Int
  , hairpinLookup :: Map Primary Int
  , hairpinGGG :: Int
  , hairpinCslope :: Int
  , hairpinCintercept :: Int
  , hairpinC3 :: Int
  , bulgeL :: PrimArray DIM1 Int
  , bulgeSingleC :: Int
  , iloop1x1 :: PrimArray PPNN  Int
  , iloop2x1 :: PrimArray PPNNN Int
  , iloop2x2 :: PrimArray PPNNNN Int
  , iloopMM :: PrimArray PNN Int
  , iloop2x3MM :: PrimArray PNN Int
  , iloop1xnMM :: PrimArray PNN Int
  , iloopL :: PrimArray DIM1 Int
  , multiMM :: PrimArray PNN Int
  , ninio :: Int
  , maxNinio :: Int
  , multiOffset :: Int
  , multiNuc :: Int
  , multiHelix :: Int
  , multiAsym :: Int
  , multiStrain :: Int
  , extMM :: PrimArray PNN Int
  , coaxial :: PrimArray PP Int -- no intervening unpaired nucleotides
  , coaxStack :: PrimArray PNN Int
  , tStackCoax :: PrimArray PNN Int
  , largeLoop :: Int
  , termAU :: Int
  , intermolecularInit :: Int
  } -- deriving (Read,Show)

