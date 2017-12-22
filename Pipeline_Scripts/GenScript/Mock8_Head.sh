./Dada2Pipe.sh ../../Config/dada2/Mock8Dada2.txt high &
./UsearchPipeF.sh ../../Config/unoise/Mock8Unoise.txt high &
./DeblurPipeF.sh ../../Config/deblur/Mock8Deblur.txt high &
wait
./Dada2Pipe.sh ../../Config/dada2/Mock8Dada2.txt med &
./UsearchPipeF.sh ../../Config/unoise/Mock8Unoise.txt med &
./DeblurPipeF.sh ../../Config/deblur/Mock8Deblur.txt med &
wait
./Dada2Pipe.sh ../../Config/dada2/Mock8Dada2.txt low &
./UsearchPipeF.sh ../../Config/unoise/Mock8Unoise.txt low &
./DeblurPipeF.sh ../../Config/deblur/Mock8Deblur.txt low &
wait
