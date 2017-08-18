for i in 0*.pnm
do
 if [ ! -e "${i}.jpg" ]
 then
  X=` echo $i | cut -f 2 -d_ --output-delimiter=\  `
  mode=` echo $i | cut -f 3 -d_ --output-delimiter=\ `
  gam=` echo $i | cut -f 4 -d_ | sed s/\.pnm//  `
  echo "${i}"
  convert -gamma 1.8 -resize 512x512 -font myriadb -fill white -pointsize 56 -draw "text 20,70 ${X}" -pointsize 28 -draw "text 380,500 ${mode}" -draw "text 20,500 ${gam}" "${i}" "${i}.jpg"
 fi
done
