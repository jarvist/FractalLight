all: glFL 

glFL: glFL.c cavity_sim.c
	gcc -o glFL glFL.c  -funroll-all-loops -funsafe-math-optimizations -O4 -finline -lm -lGL -lglut -lfftw3

mac: glFL.c cavity_sim.c
	gcc -o glFL glFL.c -funroll-all-loops -funsafe-math-optimizations -O4 -finline -lm -lfftw3 -framework GLUT -framework OpenGL 

commandlineFL: cavity_sim.c
	gcc -o commandlineFL commandlineFL.c  -funroll-all-loops -O999 -finline -lm -lfftw3

youtube-render:
	ffmpeg -framerate 10 -i VF_%08d.ppm -s:v 800x800 -c:v libx264 -profile:v high -crf 23 -pix_fmt yuv420p -r 10 VF_youtube.mp4

youtube-render-hd: # It's 2017!
	ffmpeg -framerate 30 -i pic-%06d.pbm -s:v 1080x1080 \
		-c:v libx264 -profile:v high -crf 23 -pix_fmt yuv420p -r 30 VF_youtube.mp4

youtube-render-alternative: 
	ffmpeg -framerate 10 -i pic-%06d.pbm -s:v 1080x1080 \
		-c:v libx264 -preset slow -tune animation -r 10 VF_youtube-alternative.mp4
