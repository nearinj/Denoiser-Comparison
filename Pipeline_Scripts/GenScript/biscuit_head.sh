./Dada2Pipe.sh ../../Config/dada2/BiscuitDada2.txt med &
./UsearchPipe.sh ../../Config/unoise/BiscuitUnoise.txt med &
./DeblurPipeP.sh ../../Config/deblur/BiscuitDeblur.txt med &
wait

