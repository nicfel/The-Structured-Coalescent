## Folder Description

This folder contains the 2 files esco.jar and mtt.jar. These two jar files are precompiled runnable BEAST2 jars containing the esco package (esco.jar) and the mtt + esco package (mtt.jar)


### How to run an xml

If your esco.jar file and your.xml file are in the same folder, navigate to that folder in a UNIX terminal. Then run your.xml as follows:

~~~
java -jar esco.jar your.xml
~~~

In windows command, it could work the same way but it isn't tested so far.

The jar file mtt.jar contains MultiTypeTree and esco to allow MTT to use some java classes in esco.
