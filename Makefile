CC = mpicc -fopenmp

INCDIR=src
OBJDIR=src

OBJECTS= main3.o $(OBJDIR)/io_init.o $(OBJDIR)/ordered_update.o $(OBJDIR)/static_update.o $(OBJDIR)/wave_update.o
CFLAGS = -c -I$(INCDIR) 

main3.x: $(OBJECTS)
		$(CC) $(OBJECTS) -o $@
main3.o: main3.c
		$(CC) $(CFLAGS) main3.c
$(OBJDIR)/%.o: $(INCDIR)/%.c
		$(CC) $(CFLAGS) $^ -o $@
clean:
		rm -rf $(OBJDIR)/*.o *.o
