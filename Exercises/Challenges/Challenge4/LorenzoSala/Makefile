CXX = $(mkCxxCompiler)
CXXFLAGS += -std=c++11 -Wall -pedantic

EXEC1 = main_poly
SRCS1 = main_poly.cpp polynomial.cpp

EXEC2 = main_poly_template
SRC2 = main_poly_template.cpp


all: $(EXEC1) $(EXEC2)

normal: $(EXEC1)

template: $(EXEC2)

$(EXEC1): $(SRCS1)

$(EXEC2): $(SRCS2)

clean:
	-rm *.o $(EXEC1) $(EXEC2) *.*~ *~
 
.PHONY: all clean normal template

