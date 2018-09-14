function [retval,data]=read_float_array(fid, name, ndims, dims1)
%[retval,data]=read_float_array(fid, name, ndims, dims)
% Read a float array from a dat file
%
% A Milstein 2/2001

dims=[ndims dims1];

status=find_array_with_name(fid, name, dims);
if (status ~= 0)
  retval=-1;
  return ;
end

num_elements=prod(dims1);

data=zeros(num_elements,1);

for i=1:num_elements
  [status, data_val]=read_float(fid);
  if (status ~= 0)
    retval=-1;
    return ;
  end
  data(i)=data_val;

end  

data=multitranspose(reshape(data,dims1(ndims:-1:1)));

retval=0; 

