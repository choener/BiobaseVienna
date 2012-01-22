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
import Data.PrimitiveArray.Unboxed.Zero
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
  { stack :: Arr0 PP Int
  , dangle3 :: Arr0 PN Int
  , dangle5 :: Arr0 PN Int
  , hairpinL :: Arr0 DIM1 Int
  , hairpinMM :: Arr0 PNN Int
  , hairpinLookup :: Map Primary Int
  , hairpinGGG :: Int
  , hairpinCslope :: Int
  , hairpinCintercept :: Int
  , hairpinC3 :: Int
  , bulgeL :: Arr0 DIM1 Int
  , bulgeSingleC :: Int
  , iloop1x1 :: Arr0 PPNN  Int
  , iloop2x1 :: Arr0 PPNNN Int
  , iloop2x2 :: Arr0 PPNNNN Int
  , iloopMM :: Arr0 PNN Int
  , iloop2x3MM :: Arr0 PNN Int
  , iloop1xnMM :: Arr0 PNN Int
  , iloopL :: Arr0 DIM1 Int
  , multiMM :: Arr0 PNN Int
  , ninio :: Int
  , maxNinio :: Int
  , multiOffset :: Int
  , multiNuc :: Int
  , multiHelix :: Int
  , multiAsym :: Int
  , multiStrain :: Int
  , extMM :: Arr0 PNN Int
  , coaxial :: Arr0 PP Int -- no intervening unpaired nucleotides
  , coaxStack :: Arr0 PNN Int
  , tStackCoax :: Arr0 PNN Int
  , largeLoop :: Int
  , termAU :: Int
  , intermolecularInit :: Int
  } -- deriving (Read,Show)

