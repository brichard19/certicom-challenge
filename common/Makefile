NAME=common
SRC=*.cpp
OBJS=$(SRC:.cpp=.o)

all:    ${SRC}
	for file in ${SRC} ; do\
		${CXX} -c $$file ${INCLUDE} ${CFLAGS};\
	done
	mkdir -p ${LIB_DIR}
	ar rvs ${LIB_DIR}/lib$(NAME).a *.o

clean:
	rm -rf *.o