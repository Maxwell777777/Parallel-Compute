#!/bin/sh

/opt/module/hadoop-2.8.0/bin/hadoop jar weather.jar ./London2013 ./London2013/output;
/opt/module/hadoop-2.8.0/bin/hadoop jar weather.jar ./Mumbai2013 ./Mumbai2013/output;
/opt/module/hadoop-2.8.0/bin/hadoop jar weather.jar ./NewYork2013 ./NewYork2013/output;
/opt/module/hadoop-2.8.0/bin/hadoop jar weather.jar ./SFO2012 ./SFO2012/output;
/opt/module/hadoop-2.8.0/bin/hadoop jar weather.jar ./SFO2013 ./SFO2013/output;


