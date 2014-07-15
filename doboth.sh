
source clean.sh

time root -b > an.log 2>&1 <<EOI
.L loop.C+
loop()
EOI

time root -b > an.log 2>&1 <<EOI
.L loopGMI.C+
loopGMI()
EOI
