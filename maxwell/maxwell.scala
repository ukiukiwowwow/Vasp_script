object maxwll {
   
   def main(args: Array[String]):Unit={
      import math._
      import scala.util.Random
      println("Hello, world!")
      val h=1.7*pow(10,-27)
      def makeMaxwellDistribution(M:Double,T:Int):List[Double]={
	    val pre=pow( (M/(2*Pi*1.38*pow(10,-23)*T)) , 3/2)
	    def process(v:Int):List[Double]={
	        val ev=exp(- M*pow(v,2)/(2*1.38*pow(10,-23)*T))
		    if(v==4000){
			    (4*Pi*pow(v,2)*pre*ev)::Nil
		    } else{
			    (4*Pi*pow(v,2)*pre*ev)::process(v+1)
		    }
	    }
	    process(0)
    }
      val l=makeMaxwellDistribution(h,300)
      val r=new Random()
      def applyMaxwell(list:List[Double]):Int={
          val xmin=0
          val xmax=4000
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
           case cur::cdr =>  List(cur(0)-X,cur(1)-Y,cur(2)-Z)::makeamendlist(cdr,X,Y,Z)
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
   
   
   
}
}

