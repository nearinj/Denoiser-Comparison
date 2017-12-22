./Dada2Pipe.sh ../../Config/dada2/Mock9Dada2.txt high &
./UsearchPipeF.sh ../../Config/unoise/Mock9Unoise.txt high &
./DeblurPipeF.sh ../../Config/deblur/Mock9Deblur.txt high &
wait
./Dada2Pipe.sh ../../Config/dada2/Mock9Dada2.txt med &
./UsearchPipeF.sh ../../Config/unoise/Mock9Unoise.txt med &
./DeblurPipeF.sh ../../Config/deblur/Mock9Deblur.txt med &
wait
./Dada2Pipe.sh ../../Config/dada2/Mock9Dada2.txt low &
./UsearchPipeF.sh ../../Config/unoise/Mock9Unoise.txt low &
./DeblurPipeF.sh ../../Config/deblur/Mock9Deblur.txt low &
wait
