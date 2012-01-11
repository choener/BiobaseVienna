{-# LANGUAGE TemplateHaskell #-}

-- | A set of default Turner 2004 parameters in ViennaRNA format.

module Biobase.Vienna.Default where

import Data.FileEmbed
import qualified Data.ByteString.Char8 as B

import Biobase.Vienna
import Biobase.Vienna.ImportPar



-- |

turnerRNA2004 :: (Vienna2004,Vienna2004)
turnerRNA2004 = fromByteString rnaTurner2004par

-- | embedded parameter file for Turner2004 parameters

rnaTurner2004par = $(embedFile "parfiles/rna_turner2004.par")
