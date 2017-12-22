./Dada2Pipe.sh ../../Config/dada2/ZymockDada2.txt high &
./UsearchPipe.sh ../../Config/unoise/ZymockUnoise.txt high &
./DeblurPipeP.sh ../../Config/deblur/ZymockDeblur.txt high &
wait
./Dada2Pipe.sh ../../Config/dada2/ZymockDada2.txt med &
./UsearchPipe.sh ../../Config/unoise/ZymockUnoise.txt med &
./DeblurPipeP.sh ../../Config/deblur/ZymockDeblur.txt med &
wait
./Dada2Pipe.sh ../../Config/dada2/ZymockDada2.txt low &
./UsearchPipe.sh ../../Config/unoise/ZymockUnoise.txt low &
./DeblurPipeP.sh ../../Config/deblur/ZymockDeblur.txt low &
wait
