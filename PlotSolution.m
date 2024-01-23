% Функция PlotSolution принимает решение (sol) и модель (model) в качестве входных данных
function PlotSolution(sol,model)

    % Извлечение параметров модели
    xs=model.xs;
    ys=model.ys;
    xt=model.xt;
    yt=model.yt;
    
    xobs=model.xobs;
    yobs=model.yobs;
    robs=model.robs;
    
    % Извлечение свойств пути из решения
    XS=sol.XS;
    YS=sol.YS;
    xx=sol.xx;
    yy=sol.yy;
    
    % Создание углов для отрисовки препятствий
    theta=linspace(0,2*pi,500);
    
    % Отрисовка статических препятствий
    for k=1:numel(xobs)
        fill(xobs(k)+robs(k)*cos(theta),yobs(k)+robs(k)*sin(theta),[0.5 0.7 0.8]);
        hold on;
    end
    
    % Отрисовка плана движения
    plot(xx, yy, 'k', 'LineWidth', 2);

    % Отрисовка точек управления
    plot(XS, YS, 'ro');

    % Отрисовка начальной точки (исходное положение)
    plot(xs, ys, '^', 'MarkerSize', 12, 'MarkerFaceColor', 'b');

    % Отрисовка конечной цели (конечное положение)
    plot(xt, yt, 'ko', 'MarkerSize', 16, 'MarkerFaceColor', 'g');

    % Добавление надписей
    text(xs, ys, ' Начальная точка', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
    text(xt, yt, ' Конечная цель', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');

    hold off;
    grid on;
    axis equal;


end
