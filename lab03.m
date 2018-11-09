best = -1;
for j = 1:10
    fun = @(x) x*sin(10*pi*x)+1;
    pop_size = 100;
    x = randi(2,pop_size,22, 'double') - 1;
    x_f = bi2de( x(:,:) );
    x_f = x_f(:)/(2^22)*3 - 1;
    %y = bi2de(x(1,:))
    %solution = ga(@fun, x_f);
    output = zeros(pop_size,1);
    for i = 1:pop_size
        output(i) = fun(x_f(i));
    end
    hist(output,pop_size);
    if (max(output) > best )
        best = max(output);
    end
end
best