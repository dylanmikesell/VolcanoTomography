#! /bin/sh

# 41 frequencies
snum=81

snumiter=1
while [ $snumiter -le $snum ]
do

echo $snumiter
cp runFiles/pronto$snumiter.run pronto.run
cp obsData/data$snumiter.obs data.obs
pronto.com > outFiles/pronto$snumiter.out
cp vel.final velData/vel$snumiter.final
cp ray.final rayData/ray$snumiter.final

snumiter=`bc -l <<END
             $snumiter + 1
END`
done

