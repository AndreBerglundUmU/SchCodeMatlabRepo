function plotProfile(x,fun,absCrit)
    if absCrit
        plot(x,abs(real(fun)),x,abs(imag(fun)),x,abs(fun))
        legend('|real(u(x))|','|imag(u(x))|','|u(x)|')
    else
        plot(x,real(fun),x,imag(fun),x,abs(fun))
        legend('real(u(x))','imag(u(x))','|u(x)|')
    end
end