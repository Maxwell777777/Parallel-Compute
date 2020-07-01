import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.FloatWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Mapper;
import java.io.IOException;

public class MyMapper extends Mapper<LongWritable,Text,Text,FloatWritable> {
    public void map(LongWritable key,Text value,Context context)throws IOException,InterruptedException{

        String line = value.toString();  // i.e. line = ""1","2013-01-01 00:10:00",75.2"

        String[] str = line.split(",");  // i.e str = [""1"", ""2013-01-01 00:10:00"", "75.2"]

        String[] useForKey = str[1].split(" ");  // i.e useForKey = [""2013-01-01", "00:10:00""]

        float f = Float.parseFloat(str[2]);  //.substring(0, str[2].length() - 1)

        context.write(new Text(useForKey[0]),new FloatWritable(f));  // i.e k2,v2 = ""2013-01-01",75.2
    }
}