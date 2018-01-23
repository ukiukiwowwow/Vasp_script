object maxwll {
   
   def main(args: Array[String]):Unit={
      import math._
      import java.io.PrintWriter
      import scala.util.Random
      import com.quantifind.charts.Highcharts._
      println("Hello, world!")
      val h=1.7*pow(10,-27)*2
      val Si = 28*pow(10,-27)*2
      val O = 16 * pow(10,-27)
      def makeMaxwellDistribution(M:Double,T:Int):List[Double]={
	    val pre=pow( (M/(2*Pi*1.38*pow(10,-23)*T)) , 3/2)
	    def process(v:Int):Stream[Double]={
	    	@scala.annotation.tailrec
	        val ev=exp(- M*pow(v,2)/(2*1.38*pow(10,-23)*T))
		    if(v==4001){
			    //(4*Pi*pow(v,2)*pre*ev)#::Nil
			    Stream.empty
		    } else{
			    (4*Pi*pow(v,2)*pre*ev)#::process(v+1)
		    }
	    }
	    process(0).toList
    }
      val l=makeMaxwellDistribution(Si,300)
      val c =makeMaxwellDistribution(Si,100)
      val r=new Random()
      def applyMaxwell(list:List[Double]):Int={
          val xmin=0
          val xmax=2000
          val ymax=list.max
          val xrand=r.nextDouble()
          val yrand=r.nextDouble()
          val x:Int=(xmax*xrand).toInt
          val y=ymax*yrand
          if(y<=list(x)){
          	return x
          } else{
          	applyMaxwell(list)
          }
      }
   
   
   val N = 72
   def append(n:Int):List[Int]={
       n match {
           case 0 => Nil
           case _ => applyMaxwell(l)::append(n-1)
       }
   }
   val v=append(N)
   
   def makeVelocity(V:Int):List[Double]={
       def process(ve:Double,n:Int):List[Double]={
           val vrand = r.nextDouble()
           val vsign = pow(-1,r.nextInt(2)%2)
           val v0=ve*vrand
           n match {
               case 0 => Nil
               case 1 => vsign*sqrt(ve)::process(ve-ve,n-1)
               case _ => vsign*sqrt(v0)::process(ve-v0,n-1) 
           
           }
       
       }
       
       process(sqrt(V)*sqrt(V),3)
   }
   
   def make3dlist(vlist:List[Int]):List[List[Double]]={
       vlist match{
           case x::xs => makeVelocity(x)::make3dlist(xs)
           case Nil => Nil
       }
   }
   val atomV = make3dlist(v)
   var x:Double=0
   var y:Double=0
   var z:Double=0
   for(i <- atomV){
       x+=i(0)
       y+=i(1)
       z+=i(2)
   }
   x/=N
   y/=N
   z/=N
   
   def makeamendlist(vlist:List[List[Double]],X:Double,Y:Double,Z:Double):List[List[Double]]={
       vlist match{
           case cur::cdr =>   cur.map(x => (x-X)*pow(10,-2))::makeamendlist(cdr,X,Y,Z)
           // List(cur(0)-X,cur(1)-Y,cur(2)-Z)
           case _ => Nil
       }
   }
   val flist = makeamendlist(atomV,x,y,z)
   x=0
   y=0
   z=0
   for(i <- flist){
       x+=i(0)
       y+=i(1)
       z+=i(2)
   }
   println(x,y,z)
   val file = new PrintWriter("./test.dat")
   for(i <- flist){
       println(i(0)+" "+i(1)+" "+i(2))
       file.write(i(0)+" "+i(1)+" "+i(2)+"\n")
       
   }
   file.close()
   spline((0 to 4000).toList,l)
   legend(List("300"))
   
   println(l.size,c.size)
   spline((0 to 4000).toList,c)
   legend(List("100"))
   
}
}

