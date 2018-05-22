./Dada2Pipe.sh ../../Config/dada2/Mock8.txt high &
./UsearchPipeF.sh ../../Config/unoise/Mock8Usearch.txt high &
./DeblurPipeF.sh ../../Config/deblur/Mock8.txt high &
wait
./Dada2Pipe.sh ../../Config/dada2/Mock8.txt med &
./UsearchPipeF.sh ../../Config/unoise/Mock8Usearch.txt med &
./DeblurPipeF.sh ../../Config/deblur/Mock8.txt med &
wait
./Dada2Pipe.sh ../../Config/dada2/Mock8.txt low &
./UsearchPipeF.sh ../../Config/unoise/Mock8Usearch.txt low &
./DeblurPipeF.sh ../../Config/deblur/Mock8.txt low &
wait
