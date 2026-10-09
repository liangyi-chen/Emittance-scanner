CC = g++
PKG = ibsimu-1.0.6
CXXFLAGS = -Wall -g -O0 -pthread `pkg-config --cflags $(PKG)`
LDFLAGS = -pthread -Wl,-rpath,/usr/local/lib
LDLIBS = `pkg-config --libs $(PKG)`

simu: simu.o
	$(CC) -o simu.run simu.o $(LDFLAGS) $(LDLIBS) 

simu.o: simu.cpp
	$(CC) -c -o simu.o simu.cpp $(CXXFLAGS)

plot: plot.o
	$(CC) -o plot.run plot.o $(LDFLAGS) $(LDLIBS)

plot.o: plot.cpp
	$(CC) -c -o plot.o plot.cpp $(CXXFLAGS)

save_png: save_png.o
	$(CC) -o save_png.run save_png.o $(LDFLAGS) $(LDLIBS)

save_png.o: save_png.cpp
	$(CC) -c -o save_png.o save_png.cpp $(CXXFLAGS)

clean:
	$(RM) *~ *.o simu.run plot.run save_png.run
