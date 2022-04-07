default:
	rm -fr ${FOAM_USER_LIBBIN}/lib*.so ${FOAM_USER_LIBBIN}/lib*.dylib; ./fireFoam/Allwmake

clean:
	./Allwclean
