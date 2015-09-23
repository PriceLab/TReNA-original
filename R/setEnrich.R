setEnrich = 
function( trn , set , minModSize = 15 , universe , verbose = T ) {

nTFs = ncol(trn)
if(is.null(universe)) u = rownames(trn)
if(is.character(universe)) u = intersect( universe , rownames(trn) )
tmp = matrix( NA , nrow = nTFs , ncol = 4 )
for( i in 1:nTFs ) { 
   mod = intersect( u , rownames( trn[ which( trn[,i] != 0 ) , ] ) )
   set.u = intersect( u , set )
   if( length( mod ) < minModSize ) next
   if( length( mod ) >= minModSize ) {
      t = table( u %in% mod , u %in% set.u )
      test = fisher.test( t , alternative = "greater" )
      p = test$p.value
      o  = length(intersect(mod,set.u))
      e = length(mod) * length(set.u) / length(u)
      rf = o / e
      tmp[ i , ] = c( o , e , rf , p )
   }
}

tmp = as.data.frame(tmp)
colnames(tmp) = c("overlap","expect","rf","p")
q = p.adjust(tmp$p)
tf = colnames(trn)
data.frame( tf , tmp , q )

}

