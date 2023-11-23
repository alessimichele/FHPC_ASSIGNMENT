CC = gcc -fopenmp

INCDIR=src
OBJDIR=src

OBJECTS= main.o $(OBJDIR)/io_init.o $(OBJDIR)/ordered_update.o $(OBJDIR)/static_update.o $(OBJDIR)/ordered_update_finite.o $(OBJDIR)/wave_update.o
CFLAGS = -c -I$(INCDIR) 

main.x: $(OBJECTS)
		$(CC) $(OBJECTS) -o $@
main.o: main.c
		$(CC) $(CFLAGS) main.c
$(OBJDIR)/%.o: $(INCDIR)/%.c
		$(CC) $(CFLAGS) $^ -o $@
clean:
		rm -rf $(OBJDIR)/*.o *.o
