cppFunction('NumericVector relation_vector(NumericVector x){
            int n = x.size();
            NumericVector out(0.5*n*(n-1));
            int temp=0;
            for(int i=0;i< n-1;++i){
            for(int j=i+1;j<n;++j){
            out[temp] = pow(x[i]-x[j],2);
            temp +=1;
            }
            }
            return out;
            }')

cppFunction('double vector_distance(NumericVector x, NumericVector y){
            int n=x.size();
            double out=0;
            for(int i=0;i < n; ++i){
            out += pow(x[i]-y[i],2);
            }
            out = out/n;
            return out;
            }'
)

cppFunction('NumericVector relation_matrix(NumericMatrix x){
            int nrow=x.nrow(), ncol=x.ncol();
            NumericVector out(0.5*nrow*(nrow-1));
            int temp=0;
            double result=0;
            for(int i=0;i<nrow-1;++i){
            for(int j=i+1;j<nrow;++j){
            result=0;
            for(int k=0;k < ncol; ++k){
            result += pow(x(i,k)-x(j,k),2);
            }
            result = result/ncol;
            out[temp]=result;
            temp+=1;
            }
            }
            return out;
            }')