CC=g++

all: proteins segments

proteins: proteins.o hash.o grid.o distance_metrics.o curve.o clustering.o initialization.o assignment.o update.o evaluation.o tree.o
	$(CC) -o proteins proteins.o hash.o grid.o distance_metrics.o curve.o clustering.o initialization.o assignment.o update.o evaluation.o tree.o

segments: segment.o hash.o grid.o distance_metrics.o curve.o clustering.o initialization.o assignment.o update.o evaluation.o tree.o
	$(CC) -o segments segment.o hash.o grid.o distance_metrics.o curve.o clustering.o initialization.o assignment.o update.o evaluation.o tree.o

proteins.o: proteins.cpp
	$(CC) -c proteins.cpp

segment.o: segment.cpp
	$(CC) -c segment.cpp

hash.o: hash.cpp
	$(CC) -c hash.cpp

grid.o: grid.cpp
	$(CC) -c grid.cpp

distance_metrics.o: distance_metrics.cpp
	$(CC) -c distance_metrics.cpp

curve.o: curve.cpp
	$(CC) -c curve.cpp

clustering.o: clustering.cpp
	$(CC) -c clustering.cpp

initialization.o: initialization.cpp
	$(CC) -c initialization.cpp

assignment.o: assignment.cpp
	$(CC) -c assignment.cpp

update.o: update.cpp
	$(CC) -c update.cpp

evaluation.o: evaluation.cpp
	$(CC) -c evaluation.cpp

tree.o: tree.cpp
	$(CC) -c tree.cpp


.PHONY: clean

clean:
	rm *o proteins segments

