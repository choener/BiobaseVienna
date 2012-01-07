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
  { dir :: FilePath
  , dna :: Bool
  } deriving (Show,Data,Typeable)

options = Options
  { dir = def &= args
  , dna = False &= help "apply dna prefix (default: false)"
  }

main :: IO ()
main = do
  Options{..} <- cmdArgs options
  let prefix = if dna then "dna" else ""
  trnr  <- fmap fromTurner2004 $ fromDir dir prefix "dat"
  trnrH <- fmap fromTurner2004 $ fromDir dir prefix "dh"
  putStrLn $ asPar trnr trnrH
  --mapM_ print $ int22symmetry $ iloop2x2 trnr
  return ()

{- symmetry checking
int22symmetry :: PrimArray I22K Int -> [( (I22K,Int) , (I22K,Int) )]
int22symmetry arr = filter (\(a,b) -> snd a /= snd b)
                  . map (\k -> ((k,arr!k),(swp k, arr! swp k)))
                  $ keys
  where
    keys = [ (p1,p2,(b1,b2,b3,b4))
           | p1 <- cguaP, p2 <- cguaP
           , b1 <- acgu, b2 <- acgu, b3 <- acgu, b4 <- acgu
           ]
    swp (p1,p2,(b1,b2,b3,b4)) = ( p2
                                , p1
                                , ( b3
                                  , b4
                                  , b1
                                  , b2
                                  )
                                )

type I22K = (ViennaPair,ViennaPair,(Nuc,Nuc,Nuc,Nuc))

-}
