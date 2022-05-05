.SUFFIXES:
.SUFFIXES: .o .cpp
#============================================================



TARGET1 = Conjugate_Gradien
C_OBJS1 = Conjugate_Gradient.o


ALL_SOURCES = Makefile $(C_SOURCES) $(MY_INCLUDES)


CCX = g++
CXXFLAGS = -g -Wall
#-std=c99

#============================================================
all: $(TARGET1)

.o:.cpp	$(MY_INCLUDES)
	$(CCX)  -c  $(CXXFLAGS) $<  

$(TARGET1) :   $(C_OBJS1)
	$(CCX) $(CXXFLAGS)  $^ $(LIBDIRS)  -o $@


# Implicit rules: $@ = target name, $< = first prerequisite name, $^ = name of all prerequisites
#============================================================


clean:
	rm -f $(TARGET1)   $(C_OBJS1)  *~

tar: $(ALL_SOURCES) $(DATA_FILES)
	tar cvf HW3_code.tar $(ALL_SOURCES)  $(DATA_FILES)
