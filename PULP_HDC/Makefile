PULP_OMP_APP = main
PULP_APP_OMP_SRCS = src/main.c src/associative_memory.c src/aux_functions.c
 
CORE ?= 8

PULP_CFLAGS += -DCORE=$(CORE) -O3 -g3 -Iinc
PULP_LDFLAGS = -lm
 
stackSize ?= 5120  


include $(PULP_SDK_HOME)/install/rules/pulp.mk

