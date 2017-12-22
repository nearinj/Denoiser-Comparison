./Dada2Pipe.sh ../../Config/dada2/BlueberryDada2.txt med &
./UsearchPipe.sh ../../Config/unoise/BlueberryUnoise.txt med &
./DeblurPipeP.sh ../../Config/deblur/BlueberryDeblur.txt med &
wait
