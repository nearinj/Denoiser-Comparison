./Dada2Pipe.sh ../../Config/dada2/ZymockDada2.txt high &
./UsearchPipe.sh ../../Config/unoise/ZymockUsearch.txt high &
./DeblurPipeP.sh ../../Config/deblur/Zymockdeblur.txt high &
wait
./Dada2Pipe.sh ../../Config/dada2/ZymockDada2.txt med &
./UsearchPipe.sh ../../Config/unoise/ZymockUsearch.txt med &
./DeblurPipeP.sh ../../Config/deblur/Zymockdeblur.txt med &
wait
./Dada2Pipe.sh ../../Config/dada2/ZymockDada2.txt low &
./UsearchPipe.sh ../../Config/unoise/ZymockUsearch.txt low &
./DeblurPipeP.sh ../../Config/deblur/Zymockdeblur.txt low &
wait
