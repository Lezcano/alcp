
const Fpxelem::F Fpxelem::getField()const{
    return f.lc().getField();
}

ll Fpxelem::getSize()const{
    return f.getField().getSize();
}

bool compatible(const &Fpxelem e1, const &Fpxelem e2){
    return e1.getField()==e2.getField();
}

bool operator==(const &Fpxelem lhs, ll rhs){
    return lhs.deg()==0 && lhs.lc()==lhs.getField().get(rhs);
}

bool operator==(ll lhs, const &Fpxelem rhs){
    return rhs == lhs;
}

bool operator!=(const &Fpxelem lhs, ll rhs){
    return !(lhs == rhs);
}

bool operator!=(ll lhs, const &Fpxelem rhs ){
    return !(rhs == lhs);
}

