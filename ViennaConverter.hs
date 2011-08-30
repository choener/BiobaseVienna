{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE DeriveDataTypeable #-}

-- | Convert Turner parameter files into Vienna energy files. Either into a
-- ".par" file (which is used by the Haskell RNAFold as well) or multiple ".C"
-- and ".H" files for inclusion into the ViennaRNA package.

module Main where

import System.Console.CmdArgs

import Biobase.Turner.Import
import Biobase.Vienna.Export
import Biobase.Vienna.Import



data Options = Options
  {
  } deriving (Show,Data,Typeable)

options = Options
  {
  }

main :: IO ()
main = do
  Options{..} <- cmdArgs options
  return ()
