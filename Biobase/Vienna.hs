{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE StandaloneDeriving #-}

-- | The Vienna 2004 model is a mirror of the Turner data structure. The
-- difference is that The ViennaRNA package uses an explicit annotation for RNA
-- basepairs (which this library does not anymore!)
--
-- Furthermore, all energies are stored in parts of 1/100, using Ints instead
-- of Doubles.

module Biobase.Vienna where

import qualified Data.Vector.Unboxed as VU
import qualified Data.Vector.Generic as VG
import qualified Data.Vector.Generic.Mutable as VGM
import Data.Primitive.Types

import Biobase.Turner



newtype Deka = Deka Int
  deriving (Eq,Ord,Num,Read,Show)

deriving instance Prim Deka
deriving instance VGM.MVector VU.MVector Deka
deriving instance VG.Vector VU.Vector Deka
deriving instance VU.Unbox Deka



type Vienna2004 = Turner2004Model Deka

turnerToVienna :: Turner2004 -> Vienna2004
turnerToVienna = emap (\(Energy e) -> Deka $ round $ 100 * e)

