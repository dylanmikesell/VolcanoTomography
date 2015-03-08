#! /bin/sh

# 41 frequencies
snum=62

snumiter=1
while [ $snumiter -lt $snum ]
do

cp pronto$snumiter.run pronto.run
cp data$snumiter.obs data.obs
pronto.com
cp vel.final vel$snumiter.final
cp ray.final ray$snumiter.final

snumiter=`bc -l <<END
             $snumiter + 1
END`
done

