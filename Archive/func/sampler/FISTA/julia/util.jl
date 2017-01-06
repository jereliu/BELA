# collection of link functions
function link(X, link_type)
  if link_type == "id"
    X
  elseif link_type == "pos"
    max(X, 0)
  elseif link_type == "pos2"
    max(X, 0).^2
  else
    error("illegal link function name")
  end
end

function link_deriv(X, link_type)
  if link_type == "id"
    ones(X)
  elseif link_type == "pos"
    ones(X).*(X.>0)
  elseif link_type == "pos2"
    2 * X .*(X.>0)
  else
    error("illegal link function name")
  end
end

function link_deriv2(X, link_type)
  if link_type == "id"
    zeroes(X)
  elseif link_type == "pos"
    zeroes(X)
  elseif link_type == "pos2"
    2 * ones(X) .*(X.>0)
  else
    error("illegal link function name")
  end
end
