./Dada2Pipe.sh ../../Config/dada2/Mock12.txt high &
./UsearchPipeF.sh ../../Config/unoise/Mock12Usearch.txt high &
./DeblurPipeF.sh ../../Config/deblur/Mock12.txt high &
wait
./Dada2Pipe.sh ../../Config/dada2/Mock12.txt med &
./UsearchPipeF.sh ../../Config/unoise/Mock12Usearch.txt med &
./DeblurPipeF.sh ../../Config/deblur/Mock12.txt med &
wait
./Dada2Pipe.sh ../../Config/dada2/Mock12.txt low &
./UsearchPipeF.sh ../../Config/unoise/Mock12Usearch.txt low &
./DeblurPipeF.sh ../../Config/deblur/Mock12.txt low &
wait
