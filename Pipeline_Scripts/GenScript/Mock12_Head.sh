./Dada2Pipe.sh ../../Config/dada2/Mock12Dada2.txt high &
./UsearchPipeF.sh ../../Config/unoise/Mock12Unoise.txt high &
./DeblurPipeF.sh ../../Config/deblur/Mock12Deblur.txt high &
wait
./Dada2Pipe.sh ../../Config/dada2/Mock12Dada2.txt med &
./UsearchPipeF.sh ../../Config/unoise/Mock12Unoise.txt med &
./DeblurPipeF.sh ../../Config/deblur/Mock12Deblur.txt med &
wait
./Dada2Pipe.sh ../../Config/dada2/Mock12Dada2.txt low &
./UsearchPipeF.sh ../../Config/unoise/Mock12Unoise.txt low &
./DeblurPipeF.sh ../../Config/deblur/Mock12Deblur.txt low &
wait
