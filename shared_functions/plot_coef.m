function [] = plot_coef(f, coef, lineC, lineS, width, ax, ii)

figure(ii)
semilogx(f, coef, 'Color', lineC ,'LineStyle', lineS, 'Linewidth', width);
axis(ax)
xlabel('Frequency [Hz]')
grid on
hold on

end