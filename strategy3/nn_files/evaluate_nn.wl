#!/usr/licensed/Mathematica-11.3.0/Executables/wolframscript

(* dir for nn files *)
dir = "/scratch/gpfs/dkozuch/directed_evolution/antifreeze/1hg7_evolution_nn_1/nn_files/"

(* update to avoid version probelms when importing nets *)
PacletUpdate["NeuralNetworks"]

importNets[n_] := Module[{nets, neti, netAvg, filename},
  nets = {};
  Do[{
    filename = dir <> "netAvg_p" <> ToString[i] <> ".wlnet";
    neti = Import[filename];
    AppendTo[nets, neti];
    }, {i, 1, n}];
  netAvg[x_] := Mean[#[x] & /@ nets];
  netAvg]

netAvg = importNets[5]
input = ToExpression[#]&/@{$ScriptCommandLine[[2]],$ScriptCommandLine[[3]],$ScriptCommandLine[[4]]}
score = netAvg[input]
Print[score]
