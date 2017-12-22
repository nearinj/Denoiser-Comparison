./Dada2Pipe.sh ../../Config/dada2/HMPDada2.txt high &
./UsearchPipe.sh ../../Config/unoise/HMPUnoise.txt high &
./DeblurPipeP.sh ../../Config/deblur/HMPDeblur.txt high &
wait
./Dada2Pipe.sh ../../Config/dada2/HMPDada2.txt med &
./UsearchPipe.sh ../../Config/unoise/HMPUnoise.txt med &
./DeblurPipeP.sh ../../Config/deblur/HMPDeblur.txt med &
wait
./Dada2Pipe.sh ../../Config/dada2/HMPDada2.txt low &
./UsearchPipe.sh ../../Config/unoise/HMPUnoise.txt low &
./DeblurPipeP.sh ../../Config/deblur/HMPDeblur.txt low &
wait
