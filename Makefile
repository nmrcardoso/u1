############################################################################################################
VERSION  = V1_$(shell date "+%d_%m_%Y_%T")
VERSION  = V1_$(shell date "+%d_%m_%Y_%H-%M-%S")
STANDARD = c99
############################################################################################################

INC_DIR = 
INC_LIBS = 



GCC ?= g++ -O3
# -std=c++11

all : u1




OBJS := u1.o timer.o

$(OBJDIR)%.o: $(SRCDIR)%.cpp
	$(VERBOSE)$(GCC) $(CCFLAGS) $(INC_DIR) -I. -fopenmp -MMD -MP   -c $< -o $@ 
$(OBJDIR)%.o: $(SRCDIR)%.c
	$(VERBOSE)$(GCC) $(CCFLAGS) $(INC_DIR) -I. -fopenmp -MMD -MP   -c $< -o $@ 

u1:  $(OBJS)
	$(VERBOSE)$(GCC) $(CCFLAGS)  -o $@ $+ $(INC_LIBS)  -fopenmp

clean:
	rm -f $(OBJS) u1 *.d


deps = $(OBJS:.o=.d)


pack: 
	@echo Generating Package u1_$(VERSION).tar.gz
	@tar cvfz u1_$(VERSION).tar.gz *.cpp *.h $(INCS) Makefile
	@echo Generated Package u1_$(VERSION).tar.gz

.PHONY : clean pack directories $(PROJECTNAME)

-include $(deps)
