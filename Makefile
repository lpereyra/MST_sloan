### snapshot options #######
#EXTRAS += -DLIM_VOL
#EXTRAS += -DGAL_LUM
EXTRAS += -DBRANCH_SURVIVE
EXTRAS += -DSORT_DERECHA
#EXTRAS += -DLEN_FOF_MERCHAN
#EXTRAS += -DCALCULA_MEDIA

#CC
CC     := $(OMPP) g++ $(DOMPP)
DC     := -DNTHREADS=6
CFLAGS := -Wall -O3 -fopenmp -g
GSLL   := -lgsl -lgslcblas
#VPP_INC := -I/home/lpereyra/voro++/include/voro++
#VPP_LIB := -L/home/lpereyra/voro++/lib -lvoro++
VPP_INC := -I/usr/local/include/voro++
VPP_LIB := -lvoro++
LIBS   := -lm $(GSLL) 

.PHONY : cleanall clean todo 

MAKEFILE := Makefile

OBJS := variables.o leesloan.o grid.o voronoi.o mst_kruskal.o

HEADERS := $(patsubst %.o,$.h,$(OBJS))

EXEC := mst_new.x

todo: $(EXEC)

%.o: %.c %.h $(MAKEFILE)
	$(CC) $(EXTRAS) $(VPP_INC) $(CFLAGS) $(DC) -c $<

mst_new.x: mst_new.c $(OBJS)
	$(CC) $(EXTRAS) $^ $(LIBS) $(VPP_LIB) $(CFLAGS) -o $@
	
clean:
	rm -rf $(OBJS)
	rm -rf mst_new.o
	rm -rf $(EXEC)
