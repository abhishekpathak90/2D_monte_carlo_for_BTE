function Pgon = CreateDomain(Segments)

Vertx = Segments(:,1);
Verty = Segments(:,2);

% Adding first vertex to close the loop
Vertx = [Vertx; Vertx(1,1)];
Verty = [Verty; Verty(1,1)];
P = [Vertx Verty];
Pgon = polyshape(P);