#!/bin/bash

rm -f miner-neo
echo "#!/usr/bin/Rscript" | cat - neoSourceCode.R run_neo.Rscript > miner-neo
chmod u+x miner-neo
mv miner-neo ../bin
