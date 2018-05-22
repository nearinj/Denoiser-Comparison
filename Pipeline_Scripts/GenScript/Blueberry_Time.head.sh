./Dada2Pipe.sh ../../Config/dada2/blueberry_time/5k.txt low &
./Dada2Pipe.sh ../../Config/dada2/blueberry_time/10k.txt low &
wait
./Dada2Pipe.sh ../../Config/dada2/blueberry_time/20k.txt low &
./Dada2Pipe.sh ../../Config/dada2/blueberry_time/30k.txt low &
wait
./Dada2Pipe.sh ../../Config/dada2/blueberry_time/40k.txt low &
./UsearchPipe.sh ../../Config/unoise/blueberry_time/5k.txt low &
wait
./UsearchPipe.sh ../../Config/unoise/blueberry_time/10k.txt low &
./UsearchPipe.sh ../../Config/unoise/blueberry_time/20k.txt low &
wait
./UsearchPipe.sh ../../Config/unoise/blueberry_time/30k.txt low &
./UsearchPipe.sh ../../Config/unoise/blueberry_time/40k.txt low &
wait
./DeblurPipeP.sh ../../Config/deblur/blueberry_time/5kdeblur.txt low &
./DeblurPipeP.sh ../../Config/deblur/blueberry_time/10kdeblur.txt low &
wait
./DeblurPipeP.sh ../../Config/deblur/blueberry_time/20kdeblur.txt low &
./DeblurPipeP.sh ../../Config/deblur/blueberry_time/30kdeblur.txt low &
wait
./DeblurPipeP.sh ../../Config/deblur/blueberry_time/40kdeblur.txt low &
wait
