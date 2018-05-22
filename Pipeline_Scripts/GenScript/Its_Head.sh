./Dada2Pipe.sh ../../Config/dada2/Its.txt high &
./UsearchPipeF.sh ../../Config/unoise/Mock9Usearch.txt high &
./DeblurPipeF.sh ../../Config/deblur/Mock9.txt high &
wait
./Dada2Pipe.sh ../../Config/dada2/Its.txt med &
./UsearchPipeF.sh ../../Config/unoise/Mock9Usearch.txt med &
./DeblurPipeF.sh ../../Config/deblur/Mock9.txt med &
wait
./Dada2Pipe.sh ../../Config/dada2/Its.txt low &
./UsearchPipeF.sh ../../Config/unoise/Mock9Usearch.txt low &
./DeblurPipeF.sh ../../Config/deblur/Mock9.txt low &
wait
