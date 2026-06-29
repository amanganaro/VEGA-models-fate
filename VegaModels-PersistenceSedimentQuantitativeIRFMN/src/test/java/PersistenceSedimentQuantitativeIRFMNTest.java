import insilico.core.exception.GenericFailureException;
import insilico.core.exception.InitFailureException;
import insilico.core.model.InsilicoModel;
import insilico.persistence_sediment_quantitative_irfmn.ismPersistenceSedimentQuantitativeIrfmn;
import model.ModelExecutionTest;

public class PersistenceSedimentQuantitativeIRFMNTest extends ModelExecutionTest {
    @Override
    protected InsilicoModel getModel() throws InitFailureException, GenericFailureException {
        return new ismPersistenceSedimentQuantitativeIrfmn();
    }
}
