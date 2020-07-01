import org.apache.hadoop.io.FloatWritable;
import org.apache.hadoop.mapreduce.Reducer;
import org.apache.hadoop.io.Text;
import java.io.IOException;

public class MyReduce extends Reducer<Text, FloatWritable,Text,FloatWritable> {
    public void reduce(Text key, Iterable<FloatWritable> values,Context context)throws IOException,InterruptedException{

        float max = 0;
        float min = 1024;

        for(FloatWritable value: values) {
            if(value.get()>max){
                max = value.get();
            }
            if(value.get()<min){
                min = value.get();
            }
        }
        String tem = key.toString();
        String newKey = "data is " + tem.substring(1) + " " + "max temperature is " + max + " min temperature is " + min;

        context.write(new Text(newKey), new FloatWritable(1));
    }
}