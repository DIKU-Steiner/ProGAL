

mac:
	g++ -dynamiclib -o libprimitives_mac.jnilib -I/System/Library/Frameworks/JavaVM.framework/Headers -framework JavaVM -m64 Predicates.cpp
	#g++ -fPIC -shared -c -I/System/Library/Frameworks/JavaVM.framework/Versions/Current/Headers Predicates.cpp -o libmacprimitives.so
	#/usr/lib/jvm/java-6-openjdk/include -I/usr/lib/jvm/java-6-openjdk/include/linux Predicates.cpp -o libprimitives.so

linux:
	g++ -fPIC -shared -I/usr/lib/jvm/java-6-sun/include -I/usr/lib/jvm/java-6-sun/include/linux Predicates.cpp -o libprimitives.so
